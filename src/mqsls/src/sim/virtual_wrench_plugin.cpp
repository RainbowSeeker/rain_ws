#include <gz/msgs/entity_wrench.pb.h>
#include <gz/msgs/Utility.hh>
#include <gz/transport/Node.hh>
#include <gz/sim/System.hh>
#include <gz/sim/World.hh>
#include <gz/sim/Model.hh>
#include <gz/sim/Link.hh>
#include <gz/sim/Util.hh>
#include <gz/plugin/Register.hh>
#include <gz/common/Profiler.hh>
#include <memory>
#include <string>
#include <mutex>

using namespace gz;
using namespace sim;
using namespace systems;

/// \brief Extract wrench information from a message.
/// \param[in] _ecm Entity component manager
/// \param[in] _msg Entity message. If it's a link, that link is returned. If
/// it's a model, its canonical link is returned.
/// \param[out] Force to apply.
/// \param[out] Offset of the force application point expressed in the link
/// frame.
/// \param[out] Torque to apply.
/// \param[out] Offset of the force application point expressed in the link
/// frame.
/// \return Target link entity.
Link decomposeMessage(const EntityComponentManager &_ecm,
    const msgs::EntityWrench &_msg, math::Vector3d &_force,
    math::Vector3d &_torque, math::Vector3d &_offset)
{
    if (_msg.wrench().has_force_offset())
    {
        _offset = msgs::Convert(_msg.wrench().force_offset());
    }

    if (_msg.wrench().has_force())
    {
        _force = msgs::Convert(_msg.wrench().force());
    }

    if (_msg.wrench().has_torque())
    {
        _torque = msgs::Convert(_msg.wrench().torque());
    }

    auto entity = entityFromMsg(_ecm, _msg.entity());
    if (entity == kNullEntity)
    {
        return Link();
    }

    Link link(entity);
    if (link.Valid(_ecm))
    {
        return link;
    }

    Model model(entity);
    if (model.Valid(_ecm))
    {
        return Link(model.CanonicalLink(_ecm));
    }

    gzerr << "Wrench can only be applied to a link or a model. Entity ["
            << entity << "] isn't either of them." << std::endl;
    return Link();
}

class VirtualWrenchPlugin :   
    public gz::sim::System,
    public ISystemConfigure,
    public ISystemPreUpdate
{
public:
    /// \brief Constructor
    VirtualWrenchPlugin() = default;

    /// \brief Destructor
    ~VirtualWrenchPlugin() override = default;

    // Documentation inherited
    void Configure(const Entity &_entity,
                        const std::shared_ptr<const sdf::Element> &_sdf,
                        EntityComponentManager &_ecm,
                        EventManager &/*_eventMgr*/) override
    {
        auto world = World(_entity);
        if (!world.Valid(_ecm))
        {
            gzerr << "VirtualWrenchPlugin system should be attached to a world."
                << std::endl;
            return;
        }

        // Topic to apply wrench continuously
        std::string topic = "/world/" + world.Name(_ecm).value() + "/virtual_wrench";
        if (_sdf->HasElement("topic_wrench"))
            topic = _sdf->Get<std::string>("topic_wrench");

        _node.Subscribe(topic, &VirtualWrenchPlugin::OnWrench, this);

        gzmsg << "Listening to wrench commands in [" << topic << "]"
                << std::endl;
    }

    // Documentation inherited
    void PreUpdate(const UpdateInfo &_info,
                        EntityComponentManager &_ecm) override

    {
        GZ_PROFILE("VirtualWrenchPlugin::PreUpdate");

        std::lock_guard<std::mutex> lg(_mutex);

        // Only apply wrenches when not paused
        if (_info.paused)
            return;

        // Apply wrenches
        for (const auto &[name, wrench] : _wrenches_map)
        {
            math::Vector3d force;
            math::Vector3d offset;
            math::Vector3d torque;
            auto link = decomposeMessage(_ecm, wrench, force, torque, offset);

            if (!link.Valid(_ecm))
            {
                gzerr << "Entity not found." << std::endl
                        << wrench.DebugString() << std::endl;
                continue;
            }

            link.AddWorldWrench(_ecm, force, torque, offset);

            gzmsg << "Applying wrench to [" << name << "] value: " << force
                    << std::endl;
        }
    }

private:
    void OnWrench(const msgs::EntityWrench &_msg)
    {
        std::lock_guard<std::mutex> lg(_mutex);

        if (!_msg.has_entity() || !_msg.has_wrench())
        {
            gzerr << "Missing entity or wrench in message: " << std::endl
                << _msg.DebugString() << std::endl;
            return;
        }

        _wrenches_map[_msg.entity().name()] = _msg;
    }


    /// \brief Node used to subscribe to the wrench topic
    transport::Node _node;
    std::mutex _mutex;
    std::map<std::string, msgs::EntityWrench> _wrenches_map;
};


GZ_ADD_PLUGIN(VirtualWrenchPlugin,
                    gz::sim::System,
                    VirtualWrenchPlugin::ISystemPreUpdate,
                    VirtualWrenchPlugin::ISystemConfigure)

GZ_ADD_PLUGIN_ALIAS(VirtualWrenchPlugin,
                    "virtual_wrench_plugin")